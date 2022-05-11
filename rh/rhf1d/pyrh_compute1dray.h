#include "../inputs.h"
#include "../atom.h"

typedef struct{
	int nlw, Nrays, stokes;
	double *lam, **sI, **sQ, **sU, **sV;
	double **J, **J20;
} mySpectrum;

typedef struct{
	int Nrlk;
	RLK_Line *rlk_lines;
} myRLK_Line;
 
mySpectrum rhf1d(int Ndep,
              double *rh_scale, double *rh_temp, double *rh_ne, double *rh_vz, double *rh_vmic,
              double *rh_mag, double *rh_gamma, double *rh_chi,
              double *rh_nH, int atm_scale, 
              int do_fudge, int fudge_num, double *fudge_lam, double *fudge,
              int Nloggf, int *loggf_ids, double *loggf_values);
              // double *wavetable, int Nwave); // myRLK_Line *pyrh_rlk_lines,

myRLK_Line get_RLK_lines(int argc, char *argv[]);