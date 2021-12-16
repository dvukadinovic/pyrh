#include "../inputs.h"

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, **sI, **sQ, **sU, **sV;
	double **J, **J20;
} mySpectrum;

mySpectrum rhf1d(int argc, char *argv[], int Ndep,
              double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
              double *rh_mag, double *rh_gamma, double *rh_chi,
              double **rh_nH, int atm_scale);

void get_RLK_lines(int argc, char *argv[]);