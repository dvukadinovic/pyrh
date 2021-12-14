#include "../inputs.h"

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, **sI, **sQ, **sU, **sV;
	double **J, **J20;
} mySpectrum;

mySpectrum rhf1d(int argc, char *argv[], int pyrh_Ndep,
              double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
              double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
              double **pyrh_nH, int pyrh_atm_scale);

InputData readMe(int argc, char *argv[]);