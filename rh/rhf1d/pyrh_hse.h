typedef struct{
	double **nH, *ne, *nHtot, *rho, *pg;
} myPops;

myPops hse(char* cwd, int pyrh_Ndep, double pg_top,
           double *pyrh_scale, double *pyrh_temp, double *pyrh_ne, double *pyrh_vz, double *pyrh_vmic,
           double *pyrh_mag, double *pyrh_gamma, double *pyrh_chi,
           double *pyrh_nH, double *pyrh_nHtot, int pyrh_atm_scale, 
           int do_fudge, int fudge_num, double *fudge_lam, double *fudge);

void dummy();