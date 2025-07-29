#include "../inputs.h"

typedef struct{
	double **nH, *ne, *nHtot, *rho, *pg;
} myPops;

void hse(char* cwd, int pyrh_Ndep,
           double *pyrh_scale, double *pyrh_temp, 
           double *pyrh_ne, double *pyrh_nHtot, double *pyrh_rho, double *pyrh_pg,
           int pyrh_atm_scale, 
           int fudge_num, double *fudge_lam, double *fudge,
           int Nabun, int *abundance_id, double *abundance_value);

void get_scales(char *cwd, int pyrh_Ndep,
                 double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
                 double *pyrh_nH, int pyrh_atm_scale, 
                 double lam_ref, double *tau, double *height, double *cmass,
                 int Nabun, int *abundance_id, double *abundance_value);

void get_ne_from_nH(char *cwd, 
                    int pyrh_atm_scale, int pyrh_Ndep, 
                    double *pyrh_scale, double *pyrh_temp, 
                    double *pyrh_nH, double *pyrh_ne);