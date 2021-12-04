#define bool bool_t

#include "inputs.h"

// typedef struct{
// 	int nlw, Nrays, stokes;
// 	double *lam, **sI, **sQ, **sU, **sV;
// 	double **J, **J20;
// } mySpectrum;

double** rhf1d(int argc, char *argv[], int Ndep,
              double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
              double *rh_mag, double *rh_gamma, double *rh_chi,
              double **rh_nH,
              int atm_scale);

InputData readMe(int argc, char *argv[]);