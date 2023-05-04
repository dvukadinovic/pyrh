#include "../inputs.h"

typedef struct{
	double **nH, *ne, *nHtot, *rho, *pg;
} myPops;

myPops hse(char* cwd, int pyrh_Ndep, double pg_top,
           double *pyrh_scale, double *pyrh_temp, 
           int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge);

void get_tau(char *cwd, double mu, int pyrh_Ndep, double *tau_ref,
             double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
             double *pyrh_nH, int pyrh_atm_scale, 
             double lam_ref);

void get_ne_from_nH(char *cwd, 
                    int pyrh_atm_scale, int pyrh_Ndep, 
                    double *pyrh_scale, double *pyrh_temp, 
                    double *pyrh_nH, double *pyrh_ne);