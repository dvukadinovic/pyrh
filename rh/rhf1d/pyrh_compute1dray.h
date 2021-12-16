#include "../inputs.h"
#include "../atom.h"

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, **sI, **sQ, **sU, **sV;
	double **J, **J20;
} mySpectrum;

typedef struct{
	RLK_Line *rlk_lines;
	int Nrlk;
} myRLK_Line;

mySpectrum rhf1d(int argc, char *argv[], int Ndep,
              double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
              double *rh_mag, double *rh_gamma, double *rh_chi,
              double **rh_nH, int atm_scale, int pyrh_Nrlk, RLK_Line *pyrh_rlk_lines);

int get_RLK_lines(int argc, char *argv[], RLK_Line *rlk_lines);