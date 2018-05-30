#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Parameters.h"
#include "ini.h"

void load_param() {
	f_init();
	p_init();
}

//fields initialization parameters
int type;
double eps;
double sigma;
const char* start_file;

typedef struct
{
	int type;
	double eps;
	double sigma;
	const char* file;
} init;

static int handler_fields(void* user, const char* section, const char* name, const char* value) {
	init* init_p = (init*)user;

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (MATCH("Fields Init", "type")) {
		init_p->type = atoi(value);
	}
	else if (MATCH("Fields Init", "eps")) {
		init_p->eps = atof(value);
	}
	else if (MATCH("Fields Init", "sigma")) {
		init_p->sigma = atof(value);
	}
	else if (MATCH("Fields Init", "file")) {
		init_p->file = _strdup(value);
	}
	else {
		return 0;  /* unknown section/name, error */
	}
	return 1;
}

void f_init() {
	init config;

	if (ini_parse("Config.ini", handler_fields, &config) < 0) {
		printf("Can't load 'Config.ini'\n");
		getchar();
		exit(EXIT_FAILURE);
	}
	type = config.type;
	eps = config.eps;
	sigma = config.sigma;
	start_file = config.file;
}

//Simulation parameters
//int SIZE;
double A_HORIZON;
int STEP_LIMIT;

typedef struct
{
	//int SIZE;
	double A_HORIZON;
	int STEP_LIMIT;
} s_param;

static int handler_sparam(void* user, const char* section, const char* name, const char* value) {
	s_param* s_p = (s_param*)user;

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (MATCH("Simulation Parameters", "size")) {
		//s_p->SIZE = atoi(value);
	}
	else if (MATCH("Simulation Parameters", "A_horizon")) {
		s_p->A_HORIZON = atof(value);
	}
	else if (MATCH("Simulation Parameters", "step_limit")) {
		s_p->STEP_LIMIT = atoi(value);
	}
	else {
		return 0;  /* unknown section/name, error */
	}
	return 1;
}

void p_init() {
	s_param s_par;

	if (ini_parse("Config.ini", handler_sparam, &s_par) < 0) {
		printf("Can't load 'Config.ini'\n");
		getchar();
		exit(EXIT_FAILURE);
	}
	//SIZE = s_par.SIZE;
	A_HORIZON = s_par.A_HORIZON;
	STEP_LIMIT = s_par.STEP_LIMIT;
}