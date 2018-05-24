#include "Data_collection.h"

int ncid, x_dimid, t_dimid, A_varid, delta_varid, Phi_varid, Pi_varid;
int x_varid, t_varid;
int dimids[NDIMS];
int retval;
//function to intialize the files (TEST VERSION)
void initialize_files(double *x) {
	//file creation
	if ((retval = nc_create("bigfile_test.nc", NC_CLOBBER| NC_64BIT_OFFSET, &ncid)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid, "x", SIZE, &x_dimid)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid, "time", NC_FLOAT, 1, &t_dimid, &t_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, "x", NC_FLOAT, 1, &x_dimid, &x_varid)))
		ERR(retval);
	//dimids array assignment
	dimids[0] = t_dimid;
	dimids[1] = x_dimid;

	//fields creation
	if ((retval = nc_def_var(ncid, "A", NC_FLOAT, NDIMS, dimids, &A_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, "delta", NC_FLOAT, NDIMS, dimids, &delta_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, "PHI", NC_FLOAT, NDIMS, dimids, &Phi_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, "PI", NC_FLOAT, NDIMS, dimids, &Pi_varid)))
		ERR(retval);

	//definition end
	if ((retval = nc_enddef(ncid)))
		ERR(retval);
	//Assignment of x coordinates values
	if ((retval = nc_put_var_double(ncid, x_varid, &x[0])))
		ERR(retval);
}

//function to collect data (TEST VERSION)

void data_stamp(int t_step, double t, double* A, double* delta, double* Phi, double* Pi) {
	if (t_step%TIME_STRIDE == 0) {
		size_t start[NDIMS], count[NDIMS], start_t[1], count_t[1] ;
		start[0] = t_step/TIME_STRIDE;
		start[1] = 0;
		count[0] = 1;
		count[1] = SIZE;
		start_t[0] = start[0];
		count_t[0] = count[0];
		if ((retval = nc_put_vara_double(ncid, A_varid, start, count, &A[0])))
			ERR(retval);
		if ((retval = nc_put_vara_double(ncid, delta_varid, start, count, &delta[0])))
			ERR(retval);
		if ((retval = nc_put_vara_double(ncid, Phi_varid, start, count, &Phi[0])))
			ERR(retval);
		if ((retval = nc_put_vara_double(ncid, Pi_varid, start, count, &Pi[0])))
			ERR(retval);
		if ((retval = nc_put_vara_double(ncid, t_varid, start_t, count_t, &t)))
			ERR(retval);
	}
}

//function to close the files and flush the buffers(TEST VERSION)
void close_all() {
	if ((retval = nc_close(ncid)))
		ERR(retval);
}