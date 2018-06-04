#include "Data_collection.h"
//int for error handling
int retval;

//flag for constraint writing
bool C_flag = false;

//function to inizialize all files

void initialize_files(double *x) {
	if (BIGFILE)
		initialize_bigfile(x);
	if (PROBE)
		initialize_probe(x);
	if (CONSTRAINT)
		initialize_con(x);
	if (EXTRA_PROBE)
		initialize_extra(x);
}

//bigfile data pointers
int ncid_b, x_dimid_b, t_dimid_b, A_varid_b, delta_varid_b, Phi_varid_b, Pi_varid_b;
int x_varid_b, t_varid_b;
int dimids_b[NDIMS];

//function to intialize bigfile
void initialize_bigfile(double *x) {
	//file creation
	if ((retval = nc_create("Output/ALL.nc", NC_CLOBBER | NC_64BIT_OFFSET, &ncid_b)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid_b, "time", NC_UNLIMITED, &t_dimid_b)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid_b, "x", B_SIZE, &x_dimid_b)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid_b, "time", NC_FLOAT, 1, &t_dimid_b, &t_varid_b)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_b, "x", NC_FLOAT, 1, &x_dimid_b, &x_varid_b)))
		ERR(retval);
	//dimids_b array assignment
	dimids_b[0] = t_dimid_b;
	dimids_b[1] = x_dimid_b;

	//fields creation
	if ((retval = nc_def_var(ncid_b, "A", NC_FLOAT, NDIMS, dimids_b, &A_varid_b)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_b, "delta", NC_FLOAT, NDIMS, dimids_b, &delta_varid_b)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_b, "PHI", NC_FLOAT, NDIMS, dimids_b, &Phi_varid_b)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_b, "PI", NC_FLOAT, NDIMS, dimids_b, &Pi_varid_b)))
		ERR(retval);

	//definition end
	if ((retval = nc_enddef(ncid_b)))
		ERR(retval);
	//Assignment of x coordinates values
	size_t start[1], count[1];
	ptrdiff_t stride[1], map[1];
	start[0] = 0;
	count[0] = B_SIZE;
	stride[0] = 1;
	map[0] = SPACE_STRIDE;

	if ((retval = nc_put_varm_double(ncid_b, x_varid_b, start, count, stride, map, &x[0])))
		ERR(retval);
}

//Probe files data pointers
int ncid[N_PROBE], t_dimid[N_PROBE], x_dimid[N_PROBE], A_varid[N_PROBE], delta_varid[N_PROBE], Phi_varid[N_PROBE], Pi_varid[N_PROBE];
int x_varid[N_PROBE], t_varid[N_PROBE];
int dimids[N_PROBE][NDIMS];

//funtion to initialize the probe files
void initialize_probe(double* x) {
	int i;
	int p;
	char filename[60];
	for (i = 0; i < N_PROBE; i++) {
		//file creation
		sprintf(filename, "Output/%d_%d.nc", i + 1, N_PROBE);
		if ((retval = nc_create(filename, NC_CLOBBER | NC_64BIT_OFFSET, &ncid[i])))
			ERR(retval);
		//dimensions creation
		if ((retval = nc_def_dim(ncid[i], "time", NC_UNLIMITED, &t_dimid[i])))
			ERR(retval);
		if ((retval = nc_def_dim(ncid[i], "x", 1, &x_dimid[i])))
			ERR(retval);
		//coordinate variables creation
		if ((retval = nc_def_var(ncid[i], "time", NC_FLOAT, 1, &t_dimid[i], &t_varid[i])))
			ERR(retval);
		if ((retval = nc_def_var(ncid[i], "x", NC_FLOAT, 1, &x_dimid[i], &x_varid[i])))
			ERR(retval);
		//dimids arrays assignment
		dimids[i][0] = t_dimid[i];
		dimids[i][1] = x_dimid[i];

		//fields creation
		if ((retval = nc_def_var(ncid[i], "A", NC_FLOAT, NDIMS, dimids[i], &A_varid[i])))
			ERR(retval);
		if ((retval = nc_def_var(ncid[i], "delta", NC_FLOAT, NDIMS, dimids[i], &delta_varid[i])))
			ERR(retval);
		if ((retval = nc_def_var(ncid[i], "PHI", NC_FLOAT, NDIMS, dimids[i], &Phi_varid[i])))
			ERR(retval);
		if ((retval = nc_def_var(ncid[i], "PI", NC_FLOAT, NDIMS, dimids[i], &Pi_varid[i])))
			ERR(retval);

		//definition end
		if ((retval = nc_enddef(ncid[i])))
			ERR(retval);
		//Assignment of x coordinate value
		size_t start[1];
		start[0] = 0;
		p = SIZE * i / N_PROBE;
		if ((retval = nc_put_var1_double(ncid[i], x_varid[i], start, &x[p])))
			ERR(retval);
	}
}

//Constraint file data pointers
int ncid_c, t_dimid_c, x_dimid_c, A_sq_varid_c, A_dot_varid_c, Con_varid_c;
int x_varid_c, t_varid_c;
int dimids_c[NDIMS];

//function to initialize the constraint file
void initialize_con(double *x) {
	//file creation
	if ((retval = nc_create("Output/Constraint.nc", NC_CLOBBER | NC_64BIT_OFFSET, &ncid_c)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid_c, "time", NC_UNLIMITED, &t_dimid_c)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid_c, "x", C_SIZE, &x_dimid_c)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid_c, "time", NC_FLOAT, 1, &t_dimid_c, &t_varid_c)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_c, "x", NC_FLOAT, 1, &x_dimid_c, &x_varid_c)))
		ERR(retval);
	//dimids_c array assignment
	dimids_c[0] = t_dimid_c;
	dimids_c[1] = x_dimid_c;

	//fields creation
	if ((retval = nc_def_var(ncid_c, "A^2 part", NC_FLOAT, NDIMS, dimids_c, &A_sq_varid_c)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_c, "A dot", NC_FLOAT, NDIMS, dimids_c, &A_dot_varid_c)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_c, "Constraint", NC_FLOAT, NDIMS, dimids_c, &Con_varid_c)))
		ERR(retval);
	//definition end
	if ((retval = nc_enddef(ncid_c)))
		ERR(retval);
	//Assignment of x coordinates values
	size_t start[1], count[1];
	ptrdiff_t stride[1], map[1];
	start[0] = 0;
	count[0] = C_SIZE;
	stride[0] = 1;
	map[0] = C_SPACE_STRIDE;

	if ((retval = nc_put_varm_double(ncid_c, x_varid_c, start, count, stride, map, &x[0])))
		ERR(retval);
}

//Extra probe file data pointer
int ncid_x, t_dimid_x, x_dimid_x, A_varid_x, delta_varid_x, Phi_varid_x, Pi_varid_x;
int x_varid_x, t_varid_x;
int dimids_x[NDIMS];

//function to initialize the extra probe file
void initialize_extra(double* x) {
	//file creation
	if ((retval = nc_create("Output/Extra.nc", NC_CLOBBER | NC_64BIT_OFFSET, &ncid_x)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid_x, "time", NC_UNLIMITED, &t_dimid_x)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid_x, "x", 1, &x_dimid_x)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid_x, "time", NC_FLOAT, 1, &t_dimid_x, &t_varid_x)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_x, "x", NC_FLOAT, 1, &x_dimid_x, &x_varid_x)))
		ERR(retval);
	//dimids_x array assignment
	dimids_x[0] = t_dimid_x;
	dimids_x[1] = x_dimid_x;

	//fields creation
	if ((retval = nc_def_var(ncid_x, "A", NC_FLOAT, NDIMS, dimids_x, &A_varid_x)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_x, "delta", NC_FLOAT, NDIMS, dimids_x, &delta_varid_x)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_x, "PHI", NC_FLOAT, NDIMS, dimids_x, &Phi_varid_x)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_x, "PI", NC_FLOAT, NDIMS, dimids_x, &Pi_varid_x)))
		ERR(retval);

	//definition end
	if ((retval = nc_enddef(ncid_x)))
		ERR(retval);
	//Assignment of x coordinate value
	size_t start[1];
	start[0] = 0;
	if ((retval = nc_put_var1_double(ncid_x, x_varid_x, start, &x[X_PROBE])))
		ERR(retval);
}

//function to write on all the different files
void data_stamp(int t_step, double t, double* x, double* A, double* delta, double* Phi, double* Pi, double Con[3][SIZE], double t_past) {
	if ((t_step%TIME_STRIDE == 0) && BIGFILE)
		data_stamp_bigfile(t_step, t, A, delta, Phi, Pi);
	if (t_step%PROBE_STRIDE == 0) {
		if (PROBE)
			data_stamp_probe(t_step, t, A, delta, Phi, Pi);
		if(EXTRA_PROBE)
			data_stamp_extra(t_step, t, A, delta, Phi, Pi);
	}
	if (C_flag) {
		data_stamp_con(t_step, Con, t_past);
		C_flag = false;
	}
	if ((t_step%C_TIME_STRIDE == 0) && CONSTRAINT)
		if(t_step)
			C_flag = true;
	if ((t_step%CHECKSTEP==0) && CHECKPOINT)
		data_stamp_check(t_step, t, x, A, delta, Phi, Pi);
}

//function to write data for the bigfile
void data_stamp_bigfile(int t_step, double t, double* A, double* delta, double* Phi, double* Pi) {
	size_t start[NDIMS], count[NDIMS], start_t[1];
	ptrdiff_t stride[NDIMS], map[NDIMS];
	start[0] = t_step / TIME_STRIDE;
	start[1] = 0;
	count[0] = 1;
	count[1] = B_SIZE;
	stride[0] = 1;
	stride[1] = 1;
	map[0] = 1;
	map[1] = SPACE_STRIDE;
	start_t[0] = start[0];
	if ((retval = nc_put_varm_double(ncid_b, A_varid_b, start, count, stride, map, &A[0])))
		ERR(retval);
	if ((retval = nc_put_varm_double(ncid_b, delta_varid_b, start, count, stride, map, &delta[0])))
		ERR(retval);
	if ((retval = nc_put_varm_double(ncid_b, Phi_varid_b, start, count, stride, map, &Phi[0])))
		ERR(retval);
	if ((retval = nc_put_varm_double(ncid_b, Pi_varid_b, start, count, stride, map, &Pi[0])))
		ERR(retval);
	if ((retval = nc_put_var1_double(ncid_b, t_varid_b, start_t, &t)))
		ERR(retval);
}

//function to write on the probe files
void data_stamp_probe(int t_step, double t, double* A, double* delta, double* Phi, double* Pi) {
	int i, p;
	size_t start[NDIMS], start_t[1];
	for (i = 0; i < N_PROBE; i++) {		
		start[0] = t_step / PROBE_STRIDE;
		start[1] = 0;
		start_t[0] = start[0];
		p = SIZE * i / N_PROBE;
		if ((retval = nc_put_var1_double(ncid[i], A_varid[i], start, &A[p])))
			ERR(retval);
		if ((retval = nc_put_var1_double(ncid[i], delta_varid[i], start, &delta[p])))
			ERR(retval);
		if ((retval = nc_put_var1_double(ncid[i], Phi_varid[i], start, &Phi[p])))
			ERR(retval);
		if ((retval = nc_put_var1_double(ncid[i], Pi_varid[i], start, &Pi[p])))
			ERR(retval);
		if ((retval = nc_put_var1_double(ncid[i], t_varid[i], start_t, &t)))
			ERR(retval);
	}
}

//function to write on the constraint file
void data_stamp_con(int t_step, double Con[3][SIZE], double t_past) {
	size_t start[NDIMS], count[NDIMS], start_t[1];
	ptrdiff_t stride[NDIMS], map[NDIMS];
	start[0] = ((t_step-1) / C_TIME_STRIDE)-1;
	start[1] = 0;
	count[0] = 1;
	count[1] = C_SIZE;
	stride[0] = 1;
	stride[1] = 1;
	map[0] = 1;
	map[1] = C_SPACE_STRIDE;
	start_t[0] = start[0];
	if ((retval = nc_put_varm_double(ncid_c, A_sq_varid_c, start, count, stride, map, &Con[1][0])))
		ERR(retval);
	if ((retval = nc_put_varm_double(ncid_c, A_dot_varid_c, start, count, stride, map, &Con[0][0])))
		ERR(retval);
	if ((retval = nc_put_varm_double(ncid_c, Con_varid_c, start, count, stride, map, &Con[2][0])))
		ERR(retval);	
	if ((retval = nc_put_var1_double(ncid_c, t_varid_c, start_t, &t_past)))
		ERR(retval);
}

//function to write on the extra probe file
void data_stamp_extra(int t_step, double t, double* A, double* delta, double* Phi, double* Pi) {
	size_t start[NDIMS], start_t[1];
	start[0] = t_step / PROBE_STRIDE;
	start[1] = 0;
	start_t[0] = start[0];
	if ((retval = nc_put_var1_double(ncid_x, A_varid_x, start, &A[X_PROBE])))
		ERR(retval);
	if ((retval = nc_put_var1_double(ncid_x, delta_varid_x, start, &delta[X_PROBE])))
		ERR(retval);
	if ((retval = nc_put_var1_double(ncid_x, Phi_varid_x, start, &Phi[X_PROBE])))
		ERR(retval);
	if ((retval = nc_put_var1_double(ncid_x, Pi_varid_x, start, &Pi[X_PROBE])))
		ERR(retval);
	if ((retval = nc_put_var1_double(ncid_x, t_varid_x, start_t, &t)))
		ERR(retval);
}

//data pointers and variables for checkpoint and final files
int ncid_ch, t_dimid_ch, x_dimid_ch, A_varid_ch, delta_varid_ch, Phi_varid_ch, Pi_varid_ch, size_varid_ch;
int x_varid_ch, t_varid_ch;
int dimids_ch[NDIMS];
char checkfile[60];
int checkn = 0;
int h_pos_n_id, h_pos_x_id;

//function to open, write and close a checkpoint file
void data_stamp_check(int t_step, double t, double* x, double* A, double* delta, double* Phi, double* Pi) {
	//file creation
	sprintf(checkfile, "Output/Checkpoint_%d.nc", checkn);
	checkn++;
	if ((retval = nc_create(checkfile, NC_CLOBBER | NC_64BIT_OFFSET, &ncid_ch)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid_ch, "time", 1, &t_dimid_ch)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid_ch, "x", SIZE, &x_dimid_ch)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid_ch, "time", NC_DOUBLE, 1, &t_dimid_ch, &t_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "x", NC_DOUBLE, 1, &x_dimid_ch, &x_varid_ch)))
		ERR(retval);

	//dimids_b array assignment
	dimids_ch[0] = t_dimid_ch;
	dimids_ch[1] = x_dimid_ch;

	//fields creation
	if ((retval = nc_def_var(ncid_ch, "A", NC_DOUBLE, NDIMS, dimids_ch, &A_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "delta", NC_DOUBLE, NDIMS, dimids_ch, &delta_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "PHI", NC_DOUBLE, NDIMS, dimids_ch, &Phi_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "PI", NC_DOUBLE, NDIMS, dimids_ch, &Pi_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "SIZE", NC_INT, 0, NULL, &size_varid_ch)))
		ERR(retval);

	//definition end
	if ((retval = nc_enddef(ncid_ch)))
		ERR(retval);

	//value assignment for all the variables
	if ((retval = nc_put_var_double(ncid_ch, t_varid_ch, &t)))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, x_varid_ch, &x[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, A_varid_ch, &A[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, delta_varid_ch, &delta[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, Pi_varid_ch, &Pi[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, Phi_varid_ch, &Phi[0])))
		ERR(retval);
	int sizech = SIZE;
	if ((retval = nc_put_var_int(ncid_ch, size_varid_ch, &sizech)))
		ERR(retval);

	//file closing
	if ((retval = nc_close(ncid_ch)))
		ERR(retval);
}
void data_stamp_final(double t, double* x, double* A, double* delta, double* Phi, double* Pi) {
	//file creation
	if ((retval = nc_create("Output/Final.nc", NC_CLOBBER | NC_64BIT_OFFSET, &ncid_ch)))
		ERR(retval);
	//dimensions creation
	if ((retval = nc_def_dim(ncid_ch, "time", 1, &t_dimid_ch)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid_ch, "x", SIZE, &x_dimid_ch)))
		ERR(retval);
	//coordinate variables creation
	if ((retval = nc_def_var(ncid_ch, "time", NC_DOUBLE, 1, &t_dimid_ch, &t_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "x", NC_DOUBLE, 1, &x_dimid_ch, &x_varid_ch)))
		ERR(retval);

	//dimids_b array assignment
	dimids_ch[0] = t_dimid_ch;
	dimids_ch[1] = x_dimid_ch;

	//fields creation
	if ((retval = nc_def_var(ncid_ch, "A", NC_DOUBLE, NDIMS, dimids_ch, &A_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "delta", NC_DOUBLE, NDIMS, dimids_ch, &delta_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "PHI", NC_DOUBLE, NDIMS, dimids_ch, &Phi_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "PI", NC_DOUBLE, NDIMS, dimids_ch, &Pi_varid_ch)))
		ERR(retval);
	if ((retval = nc_def_var(ncid_ch, "SIZE", NC_INT, 0, NULL, &size_varid_ch)))
		ERR(retval);

	//definition end
	if ((retval = nc_enddef(ncid_ch)))
		ERR(retval);

	//value assignment for all the variables
	if ((retval = nc_put_var_double(ncid_ch, t_varid_ch, &t)))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, x_varid_ch, &x[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, A_varid_ch, &A[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, delta_varid_ch, &delta[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, Pi_varid_ch, &Pi[0])))
		ERR(retval);
	if ((retval = nc_put_var_double(ncid_ch, Phi_varid_ch, &Phi[0])))
		ERR(retval);
	int sizech = SIZE;
	if ((retval = nc_put_var_int(ncid_ch, size_varid_ch, &sizech)))
		ERR(retval);

	//file closing
	if ((retval = nc_close(ncid_ch)))
		ERR(retval);
}

//function to close the files and flush the buffers
void close_all() {
	if (BIGFILE)
		close_bigfile();
	if (PROBE)
		close_probe();
	if (CONSTRAINT)
		close_con();
	if (EXTRA_PROBE)
		close_extra();
}

//function to close the All file
void close_bigfile() {
	if ((retval = nc_close(ncid_b)))
		ERR(retval);
}

//function to close the probe files
void close_probe() {
	int i;
	for (i = 0; i < N_PROBE; i++) {
		if ((retval = nc_close(ncid[i])))
			ERR(retval);
	}
}

//function to close the constraint file
void close_con() {
	if ((retval = nc_close(ncid_c)))
		ERR(retval);
}

//function to close 
void close_extra() {
	if ((retval = nc_close(ncid_x)))
		ERR(retval);
}

